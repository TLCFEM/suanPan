/*******************************************************************************
 * Copyright (C) 2017-2023 Theodore Chang
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

// ReSharper disable StringLiteralTypo
#include "vtkParser.h"
#include <Domain/DomainBase.h>
#include <Domain/Node.h>
#include <Element/Element.h>
#include <Toolbox/utility.h>

#ifdef SUANPAN_VTK

#include <vtkAutoInit.h>
VTK_MODULE_INIT(vtkRenderingOpenGL2)  // NOLINT(cppcoreguidelines-special-member-functions, hicpp-special-member-functions)
VTK_MODULE_INIT(vtkInteractionStyle)  // NOLINT(cppcoreguidelines-special-member-functions, hicpp-special-member-functions)
VTK_MODULE_INIT(vtkRenderingFreeType) // NOLINT(cppcoreguidelines-special-member-functions, hicpp-special-member-functions)

#include <vtkActor.h>
#include <vtkColorTransferFunction.h>
#include <vtkDataSetMapper.h>
#include <vtkDoubleArray.h>
#include <vtkLookupTable.h>
#include <vtkNamedColors.h>
#include <vtkPointData.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkScalarBarActor.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>

vtkInfo vtk_process(istringstream& command) {
    vtkInfo config;

    string keyword;

    while(!command.eof() && get_input(command, keyword))
        if(is_equal(keyword, "scale") && !get_input(command, config.scale)) config.scale = 1.;
        else if(is_equal(keyword, "deformed")) config.on_deformed = true;
        else if(is_equal(keyword, "undeformed")) config.on_deformed = false;
        else if(is_equal(keyword, "type") && get_input(command, keyword)) config.type = to_list(keyword.c_str());
        else if(is_equal(keyword, "fontsize") && !get_input(command, config.font_size)) config.font_size = 8;
        else if(is_equal(keyword, "save") && get_input(command, config.file_name)) config.save_file = true;
        else if(is_equal(keyword, "nobar")) config.colorbar = false;
        else if(is_equal(keyword, "material") && !get_input(command, config.material_type)) config.material_type = -1;
        else if(is_equal(keyword, "size")) {
            if(!get_input(command, config.canvas_size[0])) config.canvas_size[0] = 500;
            if(!get_input(command, config.canvas_size[1])) config.canvas_size[1] = 500;
        }

    return config;
}

void vtk_setup(const vtkSmartPointer<vtkUnstructuredGrid>& grid, const vtkInfo& config) {
    const auto color = vtkSmartPointer<vtkNamedColors>::New();
    const auto table = vtkSmartPointer<vtkLookupTable>::New();
    const auto func = vtkSmartPointer<vtkColorTransferFunction>::New();
    const auto mapper = vtkSmartPointer<vtkDataSetMapper>::New();
    const auto bar = vtkSmartPointer<vtkScalarBarActor>::New();
    const auto actor = vtkSmartPointer<vtkActor>::New();
    const auto interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    const auto renderer = vtkSmartPointer<vtkRenderer>::New();
    const auto window = vtkSmartPointer<vtkRenderWindow>::New();

    func->SetColorSpaceToDiverging();
    func->AddRGBPoint(0., .230, .299, .754);
    func->AddRGBPoint(1., .706, .016, .150);
    table->SetNumberOfTableValues(256);
    table->Build();
    double rgb[4] = {0., 0., 0., 1.};
    for(auto I = 0; I < 256; ++I) {
        func->GetColor(static_cast<double>(I) / 256., rgb);
        table->SetTableValue(I, rgb);
    }

    mapper->SetInputDataObject(grid);
    mapper->SetLookupTable(table);
    mapper->SetScalarRange(grid->GetPointData()->GetScalars()->GetRange());

    bar->SetLookupTable(mapper->GetLookupTable());

    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(color->GetColor3d("DodgerBlue").GetData());
    actor->GetProperty()->EdgeVisibilityOn();
    actor->GetProperty()->SetOpacity(1.);

    renderer->AddActor(actor);
    if(config.colorbar) renderer->AddActor2D(bar);
    renderer->SetBackground(color->GetColor3d("Grey").GetData());
    renderer->ResetCameraClippingRange();

    window->AddRenderer(renderer);

    interactor->SetRenderWindow(window);

    window->Render();
    window->SetSize(config.canvas_size[0], config.canvas_size[1]);
    window->SetWindowName(config.title_name.c_str());

    interactor->Start();
}

void vtk_save(vtkSmartPointer<vtkUnstructuredGrid>&& grid, const vtkInfo config) {
    // NOLINT(performance-unnecessary-value-param)
    const auto writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    writer->SetInputData(grid);
    writer->SetFileName(config.file_name.c_str());
    writer->SetFileTypeToBinary();
    writer->Write();
    suanpan_debug("Plot is written to file \"{}\".\n", config.file_name);
}

int vtk_parser(const shared_ptr<DomainBase>& domain, istringstream& command) {
    auto plot_info = vtk_process(command);
    if(!plot_info.on_deformed) plot_info.scale = 0.;

    switch(to_list(vtk_get_name(plot_info.type))) {
    case OutputType::U:
    case OutputType::V:
    case OutputType::A:
        domain->insert(make_shared<std::future<void>>(std::async(std::launch::async, vtk_plot_node_quantity, std::cref(domain), plot_info)));
        break;
    default:
        domain->insert(make_shared<std::future<void>>(std::async(std::launch::async, vtk_plot_element_quantity, std::cref(domain), plot_info)));
    }

    return SUANPAN_SUCCESS;
}

void vtk_plot_node_quantity(const shared_ptr<DomainBase>& domain, vtkInfo config) {
    auto& t_node_pool = domain->get_node_pool();
    auto& t_element_pool = domain->get_element_pool();

    config.title_name = "Plotting Nodal Quantity " + string(to_char(config.type));

    auto max_node = static_cast<unsigned>(t_node_pool.size());
    for(const auto& I : t_node_pool) max_node = std::max(max_node, I->get_tag());

    auto data = vtkSmartPointer<vtkDoubleArray>::New();
    auto node = vtkSmartPointer<vtkPoints>::New();

    data->SetNumberOfComponents(6);
    data->SetNumberOfTuples(++max_node);
    node->SetNumberOfPoints(max_node);

    for(unsigned I = 0; I < max_node; ++I) {
        node->SetPoint(I, 0., 0., 0.);
        data->SetTuple6(I, 0., 0., 0., 0., 0., 0.);
    }

    data->SetName(vtk_get_name(config.type));
    data->SetComponentName(0, "1");
    data->SetComponentName(1, "2");
    data->SetComponentName(2, "3");
    data->SetComponentName(3, "4");
    data->SetComponentName(4, "5");
    data->SetComponentName(5, "6");

    auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    grid->Allocate(static_cast<vtkIdType>(t_element_pool.size()));
    std::ranges::for_each(t_element_pool, [&](const shared_ptr<Element>& t_element) {
        t_element->SetDeformation(node, config.scale);
        t_element->GetData(data, to_list(vtk_get_name(config.type)));
        if(const auto t_cell = t_element->GetCell(); nullptr != t_cell) grid->InsertNextCell(t_cell->GetCellType(), t_cell->GetPointIds());
    });

    grid->SetPoints(node);

    if(config.store_ptr) {
        grid->GetPointData()->SetScalars(data);
        grid->GetPointData()->SetActiveScalars(vtk_get_name(config.type));
        config.grid_ptr = grid;
    }
    else if(config.save_file) {
        grid->GetPointData()->SetScalars(data);
        grid->GetPointData()->SetActiveScalars(vtk_get_name(config.type));
        auto writer = std::thread(vtk_save, std::move(grid), config);
        writer.detach();
    }
    else {
        const auto sub_data = vtkSmartPointer<vtkDoubleArray>::New();

        sub_data->SetNumberOfTuples(data->GetNumberOfTuples());
        sub_data->CopyComponent(0, data, vtk_get_index(config.type));

        grid->GetPointData()->SetScalars(sub_data);

        vtk_setup(grid, config);
    }
}

void vtk_plot_element_quantity(const shared_ptr<DomainBase>& domain, vtkInfo config) {
    auto& t_node_pool = domain->get_node_pool();
    auto& t_element_pool = domain->get_element_pool();

    config.title_name = "Plotting Element Quantity " + string(to_char(config.type));

    auto max_node = static_cast<unsigned>(t_node_pool.size());
    for(const auto& I : t_node_pool) max_node = std::max(max_node, I->get_tag());

    const auto data = vtkSmartPointer<vtkDoubleArray>::New();
    auto node = vtkSmartPointer<vtkPoints>::New();

    data->SetNumberOfComponents(6);
    data->SetNumberOfTuples(++max_node);
    node->SetNumberOfPoints(max_node);

    for(unsigned I = 0; I < max_node; ++I) {
        node->SetPoint(I, 0., 0., 0.);
        data->SetTuple6(I, 0., 0., 0., 0., 0., 0.);
    }

    data->SetName(vtk_get_name(config.type));
    data->SetComponentName(0, "1");
    data->SetComponentName(1, "2");
    data->SetComponentName(2, "3");
    data->SetComponentName(3, "4");
    data->SetComponentName(4, "5");
    data->SetComponentName(5, "6");

    mat tensor(6, max_node, fill::zeros);
    vec counter(max_node, fill::zeros);

    auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    grid->Allocate(static_cast<vtkIdType>(t_element_pool.size()));
    std::ranges::for_each(t_element_pool, [&](const shared_ptr<Element>& t_element) {
        if(-1 != config.material_type) if(const auto& t_tag = t_element->get_material_tag(); t_tag.empty() || static_cast<uword>(config.material_type) != t_tag(0)) return;
        t_element->SetDeformation(node, config.scale);
        auto& t_encoding = t_element->get_node_encoding();
        counter(t_encoding) += 1.;
        if(const auto t_data = t_element->GetData(to_list(vtk_get_name(config.type))); !t_data.empty()) tensor.cols(t_encoding) += t_data;
        if(const auto t_cell = t_element->GetCell(); nullptr != t_cell) grid->InsertNextCell(t_cell->GetCellType(), t_cell->GetPointIds());
    });

    for(unsigned I = 0; I < max_node; ++I)
        if(0. != counter(I)) {
            tensor.col(I) /= counter(I);
            data->SetTuple(I, tensor.colptr(I));
        }

    grid->SetPoints(node);

    if(config.store_ptr) {
        grid->GetPointData()->SetScalars(data);
        grid->GetPointData()->SetActiveScalars(vtk_get_name(config.type));
        config.grid_ptr = grid;
    }
    else if(config.save_file) {
        grid->GetPointData()->SetScalars(data);
        grid->GetPointData()->SetActiveScalars(vtk_get_name(config.type));
        auto writer = std::thread(vtk_save, std::move(grid), config);
        writer.detach();
    }
    else {
        const auto sub_data = vtkSmartPointer<vtkDoubleArray>::New();

        sub_data->SetNumberOfTuples(data->GetNumberOfTuples());
        sub_data->CopyComponent(0, data, vtk_get_index(config.type));

        grid->GetPointData()->SetScalars(sub_data);

        vtk_setup(grid, config);
    }
}

int vtk_get_index(const OutputType config) {
    if(config == OutputType::S11) return 0;
    if(config == OutputType::S22) return 1;
    if(config == OutputType::S33) return 2;
    if(config == OutputType::S12) return 3;
    if(config == OutputType::S23) return 4;
    if(config == OutputType::S13) return 5;
    if(config == OutputType::E11) return 0;
    if(config == OutputType::E22) return 1;
    if(config == OutputType::E33) return 2;
    if(config == OutputType::E12) return 3;
    if(config == OutputType::E23) return 4;
    if(config == OutputType::E13) return 5;
    if(config == OutputType::EE11) return 0;
    if(config == OutputType::EE22) return 1;
    if(config == OutputType::EE33) return 2;
    if(config == OutputType::EE12) return 3;
    if(config == OutputType::EE23) return 4;
    if(config == OutputType::EE13) return 5;
    if(config == OutputType::PE11) return 0;
    if(config == OutputType::PE22) return 1;
    if(config == OutputType::PE33) return 2;
    if(config == OutputType::PE12) return 3;
    if(config == OutputType::PE23) return 4;
    if(config == OutputType::PE13) return 5;
    if(config == OutputType::PE11) return 0;
    if(config == OutputType::PE22) return 1;
    if(config == OutputType::PE33) return 2;
    if(config == OutputType::PE12) return 3;
    if(config == OutputType::PE23) return 4;
    if(config == OutputType::PE13) return 5;
    if(config == OutputType::U1) return 0;
    if(config == OutputType::U2) return 1;
    if(config == OutputType::U3) return 2;
    if(config == OutputType::U4) return 3;
    if(config == OutputType::U5) return 4;
    if(config == OutputType::U6) return 5;
    if(config == OutputType::UR1) return 3;
    if(config == OutputType::UR2) return 4;
    if(config == OutputType::UR3) return 5;
    if(config == OutputType::V1) return 0;
    if(config == OutputType::V2) return 1;
    if(config == OutputType::V3) return 2;
    if(config == OutputType::V4) return 3;
    if(config == OutputType::V5) return 4;
    if(config == OutputType::V6) return 5;
    if(config == OutputType::VR1) return 3;
    if(config == OutputType::VR2) return 4;
    if(config == OutputType::VR3) return 5;
    if(config == OutputType::A1) return 0;
    if(config == OutputType::A2) return 1;
    if(config == OutputType::A3) return 2;
    if(config == OutputType::A4) return 3;
    if(config == OutputType::A5) return 4;
    if(config == OutputType::A6) return 5;
    if(config == OutputType::AR1) return 3;
    if(config == OutputType::AR2) return 4;
    if(config == OutputType::AR3) return 5;

    return 0;
}

const char* vtk_get_name(const OutputType P) {
    switch(P) {
    case OutputType::SD:
        return "SD";
    case OutputType::ED:
        return "ED";
    case OutputType::VD:
        return "VD";
    case OutputType::SS:
        return "SS";
    case OutputType::ES:
        return "ES";
    case OutputType::VS:
        return "VS";
    case OutputType::S:
    case OutputType::S11:
    case OutputType::S22:
    case OutputType::S33:
    case OutputType::S12:
    case OutputType::S23:
    case OutputType::S13:
        return "S";
    case OutputType::SINT:
        return "SINT";
    case OutputType::HYDRO:
        return "HYDRO";
    case OutputType::E:
    case OutputType::E11:
    case OutputType::E22:
    case OutputType::E33:
    case OutputType::E12:
    case OutputType::E23:
    case OutputType::E13:
        return "E";
    case OutputType::EEQ:
        return "EEQ";
    case OutputType::EINT:
        return "EINT";
    case OutputType::SP:
    case OutputType::SP1:
    case OutputType::SP2:
    case OutputType::SP3:
        return "SP";
    case OutputType::EP:
    case OutputType::EP1:
    case OutputType::EP2:
    case OutputType::EP3:
        return "EP";
    case OutputType::SINV:
        return "SINV";
    case OutputType::MISES:
        return "MISES";
    case OutputType::NMISES:
        return "NMISES";
    case OutputType::TRESC:
        return "TRESC";
    case OutputType::EE:
    case OutputType::EE11:
    case OutputType::EE22:
    case OutputType::EE33:
    case OutputType::EE12:
    case OutputType::EE23:
    case OutputType::EE13:
        return "EE";
    case OutputType::EEP:
    case OutputType::EEP1:
    case OutputType::EEP2:
    case OutputType::EEP3:
        return "EEP";
    case OutputType::EEEQ:
        return "EEEQ";
    case OutputType::PE:
    case OutputType::PE11:
    case OutputType::PE22:
    case OutputType::PE33:
    case OutputType::PE12:
    case OutputType::PE23:
    case OutputType::PE13:
        return "PE";
    case OutputType::PEP:
        return "PEP";
    case OutputType::PEEQ:
        return "PEEQ";

    case OutputType::U:
    case OutputType::UT:
    case OutputType::UR:
    case OutputType::U1:
    case OutputType::U2:
    case OutputType::U3:
    case OutputType::UR1:
    case OutputType::UR2:
    case OutputType::UR3:
    case OutputType::U4:
    case OutputType::U5:
    case OutputType::U6:
        return "U";
    case OutputType::V:
    case OutputType::VT:
    case OutputType::VR:
    case OutputType::V1:
    case OutputType::V2:
    case OutputType::V3:
    case OutputType::VR1:
    case OutputType::VR2:
    case OutputType::VR3:
    case OutputType::V4:
    case OutputType::V5:
    case OutputType::V6:
        return "V";
    case OutputType::A:
    case OutputType::AT:
    case OutputType::AR:
    case OutputType::A1:
    case OutputType::A2:
    case OutputType::A3:
    case OutputType::AR1:
    case OutputType::AR2:
    case OutputType::AR3:
    case OutputType::A4:
    case OutputType::A5:
    case OutputType::A6:
        return "A";

    case OutputType::RF:
    case OutputType::RF1:
    case OutputType::RF2:
    case OutputType::RF3:
        return "RF";
    case OutputType::DF:
    case OutputType::DF1:
    case OutputType::DF2:
    case OutputType::DF3:
    case OutputType::DF4:
    case OutputType::DF5:
    case OutputType::DF6:
        return "DF";
    case OutputType::RM:
    case OutputType::RM1:
    case OutputType::RM2:
    case OutputType::RM3:
        return "RM";
    case OutputType::RT:
        return "RT";
    case OutputType::DAMAGE:
        return "DAMAGE";
    case OutputType::DT:
        return "DT";
    case OutputType::DC:
        return "DC";
    case OutputType::KAPPAT:
        return "KAPPAT";
    case OutputType::KAPPAC:
        return "KAPPAC";
    case OutputType::VF:
        return "VF";

    case OutputType::REBARE:
        return "REBARE";
    case OutputType::REBARS:
        return "REBARS";
    case OutputType::RESULTANT:
        return "RESULTANT";
    case OutputType::AXIAL:
        return "AXIAL";
    case OutputType::SHEAR:
        return "SHEAR";
    case OutputType::MOMENT:
        return "MOMENT";
    case OutputType::TORSION:
        return "TORSION";
    case OutputType::LITR:
        return "LITR";

    default:
        return "NL";
    }
}

#endif
